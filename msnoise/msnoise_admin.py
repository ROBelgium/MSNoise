"""
MSNoise Web Administration Interface - Consolidated Version

This module provides a Flask-Admin based web interface for managing MSNoise
configurations, workflow chains, filters, and data processing jobs.

Key Features:
- Unified Configuration Management (global and workflow-specific)
- Dynamic Workflow Chain Creation and Management
- Enhanced Filter Management with Workflow Associations
- Simplified Base Classes with Modern Flask-Admin
- Generic Workflow Handling driven by Database Schema
"""

import datetime
import logging
import os
import sys
import traceback
from collections import OrderedDict

from flask import Flask, request, redirect, url_for, flash, jsonify, render_template_string
from flask_admin import Admin, BaseView, expose, AdminIndexView
from flask_admin.contrib.sqla import ModelView
from flask_admin.form import Select2Widget
from flask_admin.model import typefmt
from flask_admin import helpers as admin_helpers
from flask_babel import Babel, gettext, ngettext, lazy_gettext
from markupsafe import Markup
SafeMarkup = Markup
from sqlalchemy import and_, or_, func, desc, asc
from sqlalchemy.orm import sessionmaker
from wtforms import Form, StringField, TextAreaField, SelectField, BooleanField, IntegerField, FloatField
from wtforms.validators import DataRequired, Optional, NumberRange, ValidationError
from wtforms.widgets import TextArea, Select, CheckboxInput, NumberInput

from .api import connect, get_config, get_logger
from .api import get_config_categories_definition
from .msnoise_table_def import declare_tables


# Initialize Flask app and extensions
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'msnoise-admin-secret-key')
app.config['FLASK_ADMIN_SWATCH'] = 'cerulean'

# Babel for internationalization
babel = Babel(app)

# Database connection
db_session = connect()

# Get schema from declare_tables
schema = declare_tables()

# Logger
logger = get_logger('msnoise.admin', loglevel="INFO")


# ============================================================================
# Base Classes and Utilities
# ============================================================================

class BaseModelView(ModelView):
    """
    Enhanced base model view with common functionality for all MSNoise models.
    Simplified from the original complex implementation.
    """

    # Common settings
    can_export = True
    can_view_details = True
    page_size = 50
    can_set_page_size = True

    # Default column formatters
    column_formatters = {
        'used': lambda view, context, model, name: SafeMarkup('✓' if model.used else '✗') if hasattr(model, 'used') else '',
    }

    def __init__(self, model, session, *args, **kwargs):
        super(BaseModelView, self).__init__(model, session, *args, **kwargs)
        self.session = session

    def get_query(self):
        """Override to add common filters"""
        return self.session.query(self.model)

    def get_count_query(self):
        """Override to add common filters"""
        return self.session.query(func.count('*')).select_from(self.model)


class MSNoiseAdminIndexView(AdminIndexView):
    """Custom admin index view with MSNoise-specific dashboard"""

    @expose('/')
    def index(self):
        """Display MSNoise dashboard"""
        # Get basic statistics
        stats = {}
        try:
            stats['stations'] = db_session.query(schema.Station).filter(schema.Station.used == True).count()
            stats['filters'] = db_session.query(schema.Filter).filter(schema.Filter.used == True).count()
            stats['workflow_steps'] = db_session.query(schema.WorkflowSteps).filter(schema.WorkflowSteps.used == True).count()
            stats['jobs_total'] = db_session.query(schema.Job).count()
            stats['jobs_todo'] = db_session.query(schema.Job).filter(schema.Job.flag == 'T').count()
            stats['jobs_done'] = db_session.query(schema.Job).filter(schema.Job.flag == 'D').count()
            stats['config_params'] = db_session.query(schema.Config).filter(schema.Config.used == True).count()
        except Exception as e:
            logger.error(f"Error getting dashboard stats: {e}")
            stats = {}

        return self.render('admin/index.html', stats=stats)


# ============================================================================
# 1. Single ConfigView - Unified Configuration Management
# ============================================================================

class ConfigView(BaseModelView):
    """
    Unified configuration management for both global and workflow-specific settings.
    This replaces all the separate DvvMwcsView, DvvStretchingView, etc. classes.
    """

    column_list = ['category', 'name', 'set_number', 'value', 'description', 'used']
    column_searchable_list = ['name', 'category', 'value', 'description']
    column_filters = ['category', 'set_number', 'used']
    column_sortable_list = ['name', 'category', 'set_number']

    form_columns = ['name', 'category', 'set_number', 'value', 'param_type',
                    'default_value', 'description', 'units', 'possible_values', 'used_in', 'used']

    column_descriptions = {
        'name': 'Parameter name (e.g., maxlag, dtt_minlag)',
        'category': 'Configuration category (global, mwcs, stretching, etc.)',
        'set_number': 'Configuration set number (NULL for global configs)',
        'value': 'Current parameter value',
        'param_type': 'Data type of the parameter',
        'description': 'Description of what this parameter does',
        'used_in': 'Which processing steps use this parameter',
        'used': 'Whether this parameter is active'
    }

    form_choices = {
        'category': get_config_categories_definition(),
        'param_type': [
            ('str', 'String'),
            ('int', 'Integer'),
            ('float', 'Float'),
            ('bool', 'Boolean'),
        ]
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'value': lambda view, context, model, name: SafeMarkup(f'<code>{model.value}</code>') if model.value else '',
        'category': lambda view, context, model, name: SafeMarkup(f'<span class="label label-info">{model.category}</span>'),
        'param_type': lambda view, context, model, name: SafeMarkup(f'<span class="label label-default">{model.param_type}</span>'),
        'set_number': lambda view, context, model, name: str(model.set_number) if model.set_number is not None else 'Global',
    })

    def __init__(self, session, *args, **kwargs):
        super(ConfigView, self).__init__(schema.Config, session, *args, **kwargs)

    def create_model(self, form):
        """Create new configuration parameter with validation"""
        try:
            model = self.model()
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.add(model)
            self.session.commit()
            return model
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to create configuration: {str(e)}', 'error')
            return False

    def update_model(self, form, model):
        """Update configuration parameter with validation"""
        try:
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.commit()
            return True
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to update configuration: {str(e)}', 'error')
            return False

    @expose('/bulk_update/', methods=['GET', 'POST'])
    def bulk_update(self):
        """Bulk update configuration parameters"""
        if request.method == 'POST':
            try:
                updates = request.form.to_dict()
                for key, value in updates.items():
                    if key.startswith('config_'):
                        config_id = int(key.split('_')[1])
                        config = self.session.query(schema.Config).get(config_id)
                        if config:
                            config.value = value
                            config.updated_at = datetime.datetime.utcnow()

                self.session.commit()
                flash('Configuration updated successfully', 'success')
            except Exception as e:
                self.session.rollback()
                flash(f'Failed to update configuration: {str(e)}', 'error')

            return redirect(url_for('.index_view'))

        # GET request - show bulk update form
        configs = self.session.query(schema.Config).filter(schema.Config.used == True).all()
        return self.render('admin/config_bulk_update.html', configs=configs)


# ============================================================================
# 2. Single WorkflowStepsView - Workflow Chain Management
# ============================================================================

class WorkflowStepsView(BaseModelView):
    """
    Enhanced workflow steps management with chain validation and filter linking.
    This replaces the complex workflow-specific view classes.
    """

    column_list = ['ref', 'step_type', 'name', 'chain_order', 'parent_step', 'config_set_number', 'filters_count', 'used']
    column_searchable_list = ['step_type', 'name', 'description']
    column_filters = ['step_type', 'used', 'chain_order']
    column_sortable_list = ['ref', 'step_type', 'chain_order', 'config_set_number', 'created_at']
    form_columns = ['step_type', 'name', 'description', 'parent_step_ref', 'config_set_number', 'filters', 'used']

    list_template = "admin/model/workflow_steps_list.html"

    column_descriptions = {
        'step_type': 'Type of workflow step (mwcs, stretching, etc.)',
        'name': 'User-friendly name for this step',
        'chain_order': 'Order in the processing chain',
        'parent_step': 'Parent step in the workflow chain',
        'config_set_number': 'Configuration set number for this step',
        'filters': 'Filters associated with this workflow step (root steps only)',
        'used': 'Whether this step is active'
    }

    form_choices = {
        'step_type': get_config_categories_definition()
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'step_type': lambda view, context, model, name: SafeMarkup(f'<span class="label label-primary">{model.step_type}</span>'),
        'parent_step': lambda view, context, model, name: f'{model.parent_step.step_type} ({model.parent_step.ref})' if model.parent_step else 'Root',
        'chain_order': lambda view, context, model, name: SafeMarkup(f'<span class="badge">{model.chain_order}</span>'),
        'filters_count': lambda view, context, model, name: SafeMarkup(f'<span class="badge badge-info">{len(model.filters)}</span>') if model.filters else '0',
    })

    def __init__(self, session, *args, **kwargs):
        super(WorkflowStepsView, self).__init__(schema.WorkflowSteps, session, *args, **kwargs)

    def scaffold_form(self):
        """Override to create custom form with dynamic parent selection"""
        form_class = super(WorkflowStepsView, self).scaffold_form()

        # Replace parent_step_ref with a custom SelectField
        from wtforms import SelectField

        class CustomWorkflowStepForm(form_class):
            parent_step_ref = SelectField(
                'Parent Step',
                choices=[],
                coerce=lambda x: int(x) if x and str(x).isdigit() else None,
                validators=[Optional()]
            )

        return CustomWorkflowStepForm

    def create_form(self, obj=None):
        """Override to populate parent choices for creation"""
        form = super(WorkflowStepsView, self).create_form(obj)

        # Set initial parent choices (empty until step type is selected)
        form.parent_step_ref.choices = [('', 'No parent (root step)')]

        # Add help text for filters field
        if hasattr(form, 'filters'):
            form.filters.description = 'Only root workflow steps (without parent) can be directly linked to filters'

        return form

    def edit_form(self, obj=None):
        """Override to populate parent choices for editing"""
        form = super(WorkflowStepsView, self).edit_form(obj)

        if obj:
            # Get valid parent choices for this step type
            valid_parents = self.get_valid_parent_choices(obj.step_type, exclude_id=obj.ref)
            form.parent_step_ref.choices = valid_parents

            # Set current parent if exists
            if obj.parent_step_ref:
                form.parent_step_ref.data = obj.parent_step_ref

            # Only show filters field for root steps
            if obj.parent_step_ref is not None:
                # This is a child step, remove filters field
                if hasattr(form, 'filters'):
                    delattr(form, 'filters')
        else:
            # For new objects, start with basic choices
            form.parent_step_ref.choices = [('', 'No parent (root step)')]

        return form

    def get_valid_parent_choices(self, step_type, exclude_id=None):
        """Get valid parent choices for a given step type"""
        choices = [('', 'No parent (root step)')]

        if not step_type:
            return choices

        # Find which step types can be parents of this step type
        valid_parent_types = []
        for parent_type, child_types in schema.WorkflowSteps.VALID_CHAINS.items():
            if step_type in child_types:
                valid_parent_types.append(parent_type)

        if valid_parent_types:
            # Get workflow steps of valid parent types
            query = self.session.query(schema.WorkflowSteps).filter(
                schema.WorkflowSteps.step_type.in_(valid_parent_types),
                schema.WorkflowSteps.used == True
            )

            # Exclude current step if editing
            if exclude_id:
                query = query.filter(schema.WorkflowSteps.ref != exclude_id)

            valid_parents = query.all()

            for parent in valid_parents:
                choice_text = f'{parent.step_type} - {parent.name or "Unnamed"} (ID: {parent.ref})'
                choices.append((str(parent.ref), choice_text))

        return choices

    def create_model(self, form):
        """Create new workflow step with chain validation"""
        try:
            # Basic validation before creating model
            step_type = form.step_type.data
            parent_step_ref = form.parent_step_ref.data

            # Validate parent-child relationship
            if parent_step_ref:
                parent = self.session.query(schema.WorkflowSteps).get(parent_step_ref)
                if not parent:
                    flash('Selected parent step does not exist', 'error')
                    return False

                # Check if this is a valid chain
                valid_children = schema.WorkflowSteps.VALID_CHAINS.get(parent.step_type, [])
                if step_type not in valid_children:
                    flash(f'Invalid chain: {parent.step_type} cannot be followed by {step_type}', 'error')
                    return False

            # Create the model
            model = self.model()
            form.populate_obj(model)

            # Set chain order
            if model.parent_step_ref:
                parent = self.session.query(schema.WorkflowSteps).get(model.parent_step_ref)
                if parent:
                    model.chain_order = parent.chain_order + 1
                # Child steps should not have direct filter connections
                if hasattr(model, 'filters'):
                    model.filters = []
            else:
                model.chain_order = 0
                # Root steps must have filter connections
                if not form.filters.data:
                    flash('Root workflow steps must be connected to at least one filter', 'error')
                    return False

            self.session.add(model)
            self.session.commit()
            return model

        except Exception as e:
            self.session.rollback()
            flash(f'Failed to create workflow step: {str(e)}', 'error')
            return False

    def update_model(self, form, model):
        """Update workflow step with chain validation"""
        try:
            # Store old values
            old_parent = model.parent_step_ref
            old_step_type = model.step_type

            # Basic validation
            step_type = form.step_type.data
            parent_step_ref = form.parent_step_ref.data

            # Validate parent-child relationship
            if parent_step_ref:
                parent = self.session.query(schema.WorkflowSteps).get(parent_step_ref)
                if not parent:
                    flash('Selected parent step does not exist', 'error')
                    return False

                # Check if this is a valid chain
                valid_children = schema.WorkflowSteps.VALID_CHAINS.get(parent.step_type, [])
                if step_type not in valid_children:
                    flash(f'Invalid chain: {parent.step_type} cannot be followed by {step_type}', 'error')
                    return False

                # Check for circular reference (simple check)
                if parent_step_ref == model.ref:
                    flash('A workflow step cannot be its own parent', 'error')
                    return False

            # Update the model
            form.populate_obj(model)

            # Update chain order if parent changed
            if model.parent_step_ref != old_parent:
                if model.parent_step_ref:
                    parent = self.session.query(schema.WorkflowSteps).get(model.parent_step_ref)
                    if parent:
                        model.chain_order = parent.chain_order + 1
                else:
                    model.chain_order = 0

            # Handle filter connections
            if model.parent_step_ref is not None:
                # Child steps should not have direct filter connections
                if hasattr(model, 'filters'):
                    model.filters = []
            else:
                # Root steps should have filter connections
                if not model.filters:
                    flash('Root workflow steps should be connected to at least one filter', 'warning')

            self.session.commit()
            return True

        except Exception as e:
            self.session.rollback()
            flash(f'Failed to update workflow step: {str(e)}', 'error')
            return False

    @expose('/get_parent_options')
    def get_parent_options(self):
        """AJAX endpoint to get valid parent options for a step type"""
        step_type = request.args.get('step_type')
        current_id = request.args.get('current_id')

        if not step_type:
            return jsonify({'options': [{'value': '', 'text': 'No parent (root step)'}]})

        try:
            exclude_id = int(current_id) if current_id else None
            choices = self.get_valid_parent_choices(step_type, exclude_id=exclude_id)

            options = []
            for value, text in choices:
                options.append({'value': value, 'text': text})

            return jsonify({'options': options})

        except Exception as e:
            return jsonify({'error': str(e)})

    @expose('/chain_builder/')
    def chain_builder(self):
        """Interactive workflow chain builder"""
        filters = self.session.query(schema.Filter).filter(schema.Filter.used == True).all()
        workflow_steps = self.session.query(schema.WorkflowSteps).all()

        return self.render('admin/workflow_chain_builder.html',
                           filters=filters,
                           workflow_steps=workflow_steps,
                           valid_chains=schema.WorkflowSteps.VALID_CHAINS)

    @expose('/filter_associations/<int:step_id>')
    def filter_associations(self, step_id):
        """Show filter associations for a workflow step"""
        step = self.session.query(schema.WorkflowSteps).get_or_404(step_id)

        if step.parent_step_ref is not None:
            # Child step - show filters from root step
            root_step = step.get_root_step(self.session)
            filters = root_step.filters if root_step else []
            message = f"This child step inherits filters from root step: {root_step.name or root_step.step_type}" if root_step else "No root step found"
        else:
            # Root step - show direct filters
            filters = step.filters
            message = "This root step has the following direct filter associations:"

        return self.render('admin/workflow_step_filters.html',
                           step=step,
                           filters=filters,
                           message=message)

# ============================================================================
# 3. Enhanced FilterView - Filter Management with Workflow Associations
# ============================================================================

class FilterView(BaseModelView):
    """
    Enhanced filter management with workflow associations and validation.
    """

    column_list = ['ref', 'freqmin', 'freqmax', 'CC', 'SC', 'AC', 'used', 'workflow_steps_count']
    column_searchable_list = ['freqmin', 'freqmax']
    column_filters = ['CC', 'SC', 'AC', 'used']
    column_sortable_list = ['ref', 'freqmin', 'freqmax', 'used']

    form_columns = ['freqmin', 'freqmax', 'CC', 'SC', 'AC', 'workflow_steps', 'used']

    column_descriptions = {
        'freqmin': 'Lower frequency bound (Hz)',
        'freqmax': 'Upper frequency bound (Hz)',
        'CC': 'Compute cross-correlation between different stations',
        'SC': 'Compute cross-correlation between different components of same station',
        'AC': 'Compute auto-correlation from single station components',
        'workflow_steps': 'Root workflow steps that use this filter',
        'used': 'Whether this filter is active'
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'freqmin': lambda view, context, model, name: f'{model.freqmin:.3f}' if model.freqmin else '',
        'freqmax': lambda view, context, model, name: f'{model.freqmax:.3f}' if model.freqmax else '',
        'workflow_steps_count': lambda view, context, model, name: SafeMarkup(f'<span class="badge badge-success">{len(model.workflow_steps)}</span>') if model.workflow_steps else '0',
    })

    def __init__(self, session, *args, **kwargs):
        super(FilterView, self).__init__(schema.Filter, session, *args, **kwargs)

    def get_form(self):
        """Override to customize the workflow_steps field"""
        form = super(FilterView, self).get_form()

        # Customize the workflow_steps field to only show root steps
        if hasattr(form, 'workflow_steps'):
            form.workflow_steps.query_factory = lambda: self.session.query(schema.WorkflowSteps).filter(
                schema.WorkflowSteps.parent_step_ref.is_(None),
                schema.WorkflowSteps.used == True
            ).all()

            form.workflow_steps.render_kw = {
                'data-toggle': 'tooltip',
                'data-placement': 'top',
                'title': 'Only root workflow steps (without parent) can be directly linked to filters'
            }

        return form

    def create_model(self, form):
        """Create new filter with frequency validation"""
        try:
            model = self.model()
            form.populate_obj(model)

            # Validate frequencies
            model.validate_frequencies()

            # Validate that only root workflow steps are associated
            for step in model.workflow_steps:
                if step.parent_step_ref is not None:
                    raise ValueError(f"Cannot associate filter with child workflow step: {step.name}")

            self.session.add(model)
            self.session.commit()
            return model
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to create filter: {str(e)}', 'error')
            return False

    def update_model(self, form, model):
        """Update filter with frequency validation"""
        try:
            form.populate_obj(model)

            # Validate frequencies
            model.validate_frequencies()

            # Validate that only root workflow steps are associated
            for step in model.workflow_steps:
                if step.parent_step_ref is not None:
                    raise ValueError(f"Cannot associate filter with child workflow step: {step.name}")

            self.session.commit()
            return True
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to update filter: {str(e)}', 'error')
            return False

    @expose('/workflow_associations/<int:filter_id>')
    def workflow_associations(self, filter_id):
        """Show workflow associations for a filter"""
        filter_obj = self.session.query(schema.Filter).get_or_404(filter_id)
        workflows = filter_obj.get_processing_workflows(self.session)

        return self.render('admin/filter_workflows.html',
                           filter=filter_obj,
                           workflows=workflows)

# ============================================================================
# 4. Simplified Base Classes - Standard Model Views
# ============================================================================

class StationView(BaseModelView):
    """Simplified station management"""

    column_list = ['ref', 'net', 'sta', 'X', 'Y', 'altitude', 'coordinates', 'used']
    column_searchable_list = ['net', 'sta']
    column_filters = ['net', 'coordinates', 'used']
    column_sortable_list = ['ref', 'net', 'sta', 'X', 'Y', 'used']

    form_columns = ['net', 'sta', 'X', 'Y', 'altitude', 'coordinates',
                    'used_location_codes', 'used_channel_names', 'used']

    def __init__(self, session, *args, **kwargs):
        super(StationView, self).__init__(schema.Station, session, *args, **kwargs)


class JobView(BaseModelView):
    """Simplified job management"""

    column_list = ['ref', 'day', 'pair', 'jobtype', 'flag', 'lastmod']
    column_searchable_list = ['day', 'pair', 'jobtype']
    column_filters = ['jobtype', 'flag']
    column_sortable_list = ['ref', 'day', 'pair', 'jobtype', 'flag', 'lastmod']

    form_columns = ['day', 'pair', 'jobtype', 'flag']

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'flag': lambda view, context, model, name: SafeMarkup(
            f'<span class="label label-{"success" if model.flag == "D" else "warning" if model.flag == "I" else "default"}">'
            f'{"Done" if model.flag == "D" else "In Progress" if model.flag == "I" else "Todo"}</span>'
        ),
        'lastmod': lambda view, context, model, name: model.lastmod.strftime('%Y-%m-%d %H:%M:%S') if model.lastmod else '',
    })

    def __init__(self, session, *args, **kwargs):
        super(JobView, self).__init__(schema.Job, session, *args, **kwargs)


class DataAvailabilityView(BaseModelView):
    """Simplified data availability management"""

    column_list = ['ref', 'net', 'sta', 'loc', 'chan', 'starttime', 'endtime', 'samplerate', 'flag']
    column_searchable_list = ['net', 'sta', 'loc', 'chan', 'file']
    column_filters = ['net', 'sta', 'chan', 'flag']
    column_sortable_list = ['ref', 'net', 'sta', 'starttime', 'endtime', 'samplerate']

    # Read-only view
    can_create = False
    can_edit = False
    can_delete = False

    def __init__(self, session, *args, **kwargs):
        super(DataAvailabilityView, self).__init__(schema.DataAvailability, session, *args, **kwargs)

from flask_admin import BaseView, expose
from flask_admin.form import BaseForm
from flask import request, redirect, url_for, flash
from wtforms import Form, StringField, FloatField, IntegerField, BooleanField, SelectField, TextAreaField
from wtforms.validators import DataRequired, Optional
import datetime

class ConfigSetView(BaseView):
    """
    View for managing configuration sets as single forms.
    Allows editing all parameters in a set together.
    """

    def __init__(self, session, *args, **kwargs):
        super(ConfigSetView, self).__init__(*args, **kwargs)
        self.session = session

    @expose('/')
    def index(self):
        """Display all configuration sets grouped by category"""

        """Display all configuration sets grouped by category"""
        from .api import get_config_sets_organized

        ordered_categories = get_config_sets_organized(self.session)

        return self.render('admin/config_sets_index.html',
                           categories=ordered_categories)


    @expose('/create/<category>/', methods=['GET', 'POST'])
    def create_set(self, category):
        """Create a new configuration set for a category"""

        from .api import create_config_set
        set_number = create_config_set(self.session, category)


        if set_number is not None:
            self.session.commit()
            flash(f'Configuration set "{category}" #{set_number} created successfully', 'success')
            return redirect(url_for('.edit_set', category=category, set_number=set_number))

        return

    @expose('/view/<category>/<int:set_number>/')
    def view_set(self, category, set_number):
        """View configuration set in read-only mode"""

        # Get all configs for this category and set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Group configs by logical sections if needed
        config_groups = {}
        for config in configs:
            config.description = Markup(config.description)
            # You can add logic here to group configs by prefix or type
            group_name = "Configuration Parameters"
            if group_name not in config_groups:
                config_groups[group_name] = []

            config_groups[group_name].append(config)

        return self.render('admin/config_set_view.html',
                           category=category,
                           set_number=set_number,
                           configs=configs,
                           config_groups=config_groups)


    @expose('/edit/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def edit_set(self, category, set_number):
        """Edit a complete configuration set as a single form"""

        # Get all config parameters for this set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Create dynamic form class
        form_class = self._create_dynamic_form(configs)
        # Initialize form with request data for POST requests
        if request.method == 'POST':
            form = form_class(request.form)
        else:
            form = form_class()

        if request.method == 'POST' and form.validate():
            try:
                # Update all config values
                for config in configs:
                    field_name = f"param_{config.ref}"
                    if hasattr(form, field_name):
                        field = getattr(form, field_name)
                        new_value = field.data

                        # Handle different field types and preserve existing values if new value is empty
                        if config.param_type == 'bool':
                            # For boolean fields
                            if isinstance(field, BooleanField):
                                # WTForms BooleanField returns True/False
                                new_value = 'Y' if new_value else 'N'
                            elif isinstance(field, SelectField):
                                # SelectField for boolean with Y/N choices
                                new_value = str(new_value) if new_value else config.value
                            else:
                                new_value = str(new_value) if new_value else config.value
                        elif config.param_type == 'int':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(int(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        elif config.param_type == 'float':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(float(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        else:  # string type
                            if new_value is not None:
                                new_value = str(new_value).strip()
                                # Allow empty strings for string fields
                            else:
                                new_value = config.value

                        # Only update if we have a valid new value
                        if new_value is not None:
                            config.value = new_value
                            config.updated_at = datetime.datetime.utcnow()

                            # Validate against possible values (only for non-empty values)
                            if new_value.strip():
                                try:
                                    config.validate_value()
                                except ValueError as ve:
                                    flash(f'Validation error for {config.name}: {str(ve)}', 'error')
                                    raise ve
                            # flash(f'Updated parameter "{config.name}"="{new_value}"', 'success')
                            self.session.commit()
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} updated successfully', 'success')
                return redirect(url_for('.view_set', category=category, set_number=set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error updating configuration: {str(e)}', 'error')

        else:
            # Populate form with current values (GET request or validation failed)
            for config in configs:
                config.description = Markup(config.description)
                field_name = f"param_{config.ref}"
                if hasattr(form, field_name):
                    field = getattr(form, field_name)
                    current_value = config.value if config.value is not None else ''
                    field.description = config.description
                    if config.param_type == 'bool':
                        # Handle boolean conversion for different formats
                        if isinstance(field, BooleanField):
                            field.data = current_value.upper() in ('Y', 'TRUE', 'YES', '1', 'ON')
                        else:  # SelectField for boolean
                            field.data = current_value
                    elif config.param_type in ['int', 'float']:
                        try:
                            if current_value and str(current_value).strip():
                                if config.param_type == 'float':
                                    field.data = float(current_value)
                                else:
                                    field.data = int(current_value)
                            else:
                                field.data = None
                        except (ValueError, TypeError):
                            # If conversion fails, set as string for display
                            field.data = current_value
                    else:
                        field.data = current_value

            # Show validation errors if any
            if request.method == 'POST':
                for field_name, errors in form.errors.items():
                    for error in errors:
                        flash(f'Validation error in {field_name}: {error}', 'error')

        return self.render('admin/config_set_edit.html',
                           form=form,
                           configs=configs,
                           category=category,
                           set_number=set_number)

    @expose('/copy/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def copy_set(self, category, set_number):
        """Copy a configuration set to create a new one"""

        if request.method == 'POST':
            new_category = request.form.get('new_category', category)
            new_set_number = int(request.form.get('new_set_number', set_number + 1))

            try:
                # Check if target already exists
                existing = self.session.query(schema.Config).filter(
                    schema.Config.category == new_category,
                    schema.Config.set_number == new_set_number
                ).first()

                if existing:
                    flash(f'Configuration set "{new_category}" #{new_set_number} already exists', 'error')
                    return redirect(url_for('.index'))

                # Get source configs
                source_configs = self.session.query(schema.Config).filter(
                    schema.Config.category == category,
                    schema.Config.set_number == set_number
                ).all()

                if not source_configs:
                    flash(f'Source configuration set not found', 'error')
                    return redirect(url_for('.index'))

                # Copy configs
                for config in source_configs:
                    new_config = schema.Config(
                        name=config.name,
                        category=new_category,
                        set_number=new_set_number,
                        value=config.value,
                        param_type=config.param_type,
                        default_value=config.default_value,
                        description=config.description,
                        units=config.units,
                        possible_values=config.possible_values,
                        used_in=config.used_in,
                        used=config.used
                    )
                    self.session.add(new_config)

                self.session.commit()
                flash(f'Configuration set copied to "{new_category}" #{new_set_number}', 'success')
                return redirect(url_for('.edit_set', category=new_category, set_number=new_set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error copying configuration: {str(e)}', 'error')

        return self.render('admin/config_set_copy.html',
                           category=category,
                           set_number=set_number)

    @expose('/delete/<category>/<int:set_number>/', methods=['POST'])
    def delete_set(self, category, set_number):
        """Delete a configuration set"""

        try:
            deleted_count = self.session.query(schema.Config).filter(
                schema.Config.category == category,
                schema.Config.set_number == set_number
            ).delete()
            self.session.commit()
            if deleted_count > 0:
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} deleted ({deleted_count} parameters)', 'success')
            else:
                flash(f'Configuration set "{category}" #{set_number} not found', 'error')

        except Exception as e:
            self.session.rollback()
            flash(f'Error deleting configuration: {str(e)}', 'error')

        return redirect(url_for('.index'))

    def _create_dynamic_form(self, configs):
        """Create a dynamic WTForms class based on config parameters"""

        class DynamicConfigForm(Form):
            pass

        for config in configs:
            field_name = f"param_{config.ref}"

            # Determine field type based on param_type
            if config.param_type == 'int':
                field = IntegerField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'float':
                field = FloatField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'bool':
                # Check if we should use a dropdown or checkbox
                field = BooleanField(
                    label=config.name,
                    description=config.description
                )
            else:  # string type
                # Check if there are possible values (dropdown)
                if config.possible_values:
                    # Parse possible values and ensure current value is included
                    possible_vals = [v.strip() for v in config.possible_values.split('/') if v.strip()]

                    # Make sure the current value is in the choices (for validation)
                    current_value = config.value.strip() if config.value else ''
                    if current_value and current_value not in possible_vals:
                        possible_vals.append(current_value)

                    # Create choices as (value, label) tuples
                    choices = [(v, v) for v in possible_vals]

                    # Add empty choice for optional fields
                    if not current_value or 'optional' in (config.description or '').lower():
                        choices.insert(0, ('', '-- Select --'))

                    field = SelectField(
                        label=config.name,
                        description=config.description,
                        choices=choices,
                        validators=[Optional()]
                    )
                else:
                    # Long text or regular string
                    if config.description and len(config.description) > 100:
                        field = TextAreaField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )
                    else:
                        field = StringField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )

            # # Add units to description if available
            # if config.units:
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} (Units: {config.units})".strip()
            #
            # # Add possible values to description if available and not already a select field
            # if config.possible_values and not isinstance(field, SelectField):
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} Possible values: {config.possible_values}".strip()

            setattr(DynamicConfigForm, field_name, field)

        return DynamicConfigForm


# ============================================================================
# 5. Generic Workflow Handling - API Endpoints
# ============================================================================

class WorkflowAPIView(BaseView):
    """
    Generic workflow handling driven by database schema.
    Provides REST-like API endpoints for workflow management.
    """

    @expose('/api/workflow_types')
    def get_workflow_types(self):
        """Get available workflow step types"""
        return jsonify({
            'workflow_types': list(schema.WorkflowSteps.VALID_CHAINS.keys()),
            'valid_chains': schema.WorkflowSteps.VALID_CHAINS
        })

    @expose('/api/config_categories')
    def get_config_categories(self):
        """Get available configuration categories"""
        categories = db_session.query(schema.Config.category).all()
        return jsonify({
            'categories': [cat[0] for cat in categories]
        })

    @expose('/api/validate_workflow_chain', methods=['POST'])
    def validate_workflow_chain_api(self):
        """Validate a workflow chain via API"""
        try:
            chain_data = request.json
            steps = chain_data.get('steps', [])

            # Validate each step in the chain
            for i, step in enumerate(steps):
                if i > 0:
                    parent_type = steps[i-1]['type']
                    current_type = step['type']

                    if current_type not in schema.WorkflowSteps.VALID_CHAINS.get(parent_type, []):
                        return jsonify({
                            'valid': False,
                            'message': f'Invalid chain: {parent_type} cannot be followed by {current_type}'
                        })

            return jsonify({'valid': True, 'message': 'Chain is valid'})
        except Exception as e:
            return jsonify({'valid': False, 'message': str(e)})

    @expose('/api/create_workflow_chain', methods=['POST'])
    def create_workflow_chain_api(self):
        """Create a complete workflow chain via API"""
        try:
            chain_data = request.json
            steps = chain_data.get('steps', [])
            filter_ids = chain_data.get('filter_ids', [])

            created_steps = []
            parent_ref = None

            for i, step_data in enumerate(steps):
                step = schema.WorkflowSteps(
                    step_type=step_data['type'],
                    name=step_data.get('name', f"{step_data['type']} step {i+1}"),
                    description=step_data.get('description', ''),
                    parent_step_ref=parent_ref,
                    config_set_number=step_data.get('config_set_number', 1),
                    chain_order=i
                )

                db_session.add(step)
                db_session.flush()  # Get the ID

                # Connect filters to root step
                if i == 0 and filter_ids:
                    filters = db_session.query(schema.Filter).filter(schema.Filter.ref.in_(filter_ids)).all()
                    step.filters = filters

                created_steps.append(step)
                parent_ref = step.ref

            db_session.commit()

            return jsonify({
                'success': True,
                'message': f'Created workflow chain with {len(created_steps)} steps',
                'step_ids': [step.ref for step in created_steps]
            })
        except Exception as e:
            db_session.rollback()
            return jsonify({'success': False, 'message': str(e)})


# ============================================================================
# Custom Templates
# ============================================================================

# Custom templates for enhanced UI
WORKFLOW_CHAIN_BUILDER_TEMPLATE = """
{% extends 'admin/base.html' %}

{% block body %}
<div class="container-fluid">
    <h2>Workflow Chain Builder</h2>
    
    <div class="row">
        <div class="col-md-4">
            <div class="panel panel-default">
                <div class="panel-heading">Available Filters</div>
                <div class="panel-body">
                    {% for filter in filters %}
                    <div class="checkbox">
                        <label>
                            <input type="checkbox" name="filters" value="{{ filter.ref }}">
                            {{ filter.freqmin }}-{{ filter.freqmax }} Hz
                        </label>
                    </div>
                    {% endfor %}
                </div>
            </div>
        </div>
        
        <div class="col-md-8">
            <div class="panel panel-default">
                <div class="panel-heading">Workflow Chain</div>
                <div class="panel-body">
                    <div id="workflow-chain">
                        <!-- Chain builder interface will be here -->
                    </div>
                    <button class="btn btn-primary" onclick="createWorkflowChain()">Create Chain</button>
                </div>
            </div>
        </div>
    </div>
</div>

<script>
var validChains = {{ valid_chains|tojson }};

function createWorkflowChain() {
    // Implementation for creating workflow chain
    console.log('Creating workflow chain...');
}
</script>
{% endblock %}
"""

CONFIG_BULK_UPDATE_TEMPLATE = """
{% extends 'admin/base.html' %}

{% block body %}
<div class="container-fluid">
    <h2>Bulk Configuration Update</h2>
    
    <form method="post">
        {% for config in configs %}
        <div class="form-group">
            <label>{{ config.name }} ({{ config.category }})</label>
            <input type="text" name="config_{{ config.ref }}" value="{{ config.value }}" class="form-control">
            <small class="help-block">{{ config.description }}</small>
        </div>
        {% endfor %}
        
        <button type="submit" class="btn btn-primary">Update All</button>
    </form>
</div>
{% endblock %}
"""

MSNOISE_INDEX_TEMPLATE = """
{% extends 'admin/index.html' %}

{% block body %}
<div class="container-fluid">
    <h2>MSNoise Dashboard</h2>
    
    <div class="row">
        <div class="col-md-3">
            <div class="panel panel-primary">
                <div class="panel-heading">Stations</div>
                <div class="panel-body">
                    <h3>{{ stats.get('stations', 0) }}</h3>
                    <p>Active stations</p>
                </div>
            </div>
        </div>
        
        <div class="col-md-3">
            <div class="panel panel-info">
                <div class="panel-heading">Filters</div>
                <div class="panel-body">
                    <h3>{{ stats.get('filters', 0) }}</h3>
                    <p>Active filters</p>
                </div>
            </div>
        </div>
        
        <div class="col-md-3">
            <div class="panel panel-success">
                <div class="panel-heading">Workflows</div>
                <div class="panel-body">
                    <h3>{{ stats.get('workflow_steps', 0) }}</h3>
                    <p>Active workflow steps</p>
                </div>
            </div>
        </div>
        
        <div class="col-md-3">
            <div class="panel panel-warning">
                <div class="panel-heading">Jobs</div>
                <div class="panel-body">
                    <h3>{{ stats.get('jobs_todo', 0) }} / {{ stats.get('jobs_total', 0) }}</h3>
                    <p>Todo / Total jobs</p>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
"""

import csv
import os
from flask import request, redirect, url_for, flash, render_template_string, render_template
from flask_admin import BaseView, expose
from wtforms import Form, SelectField, StringField, TextAreaField, IntegerField, validators

# ============================================================================
# Main Application Setup
# ============================================================================
 
def create_admin_app():
    """Create and configure the Flask-Admin application"""

    # Initialize Flask-Admin
    admin = Admin(
        app,
        name='MSNoise Admin',
        template_mode='bootstrap3',
        index_view=MSNoiseAdminIndexView(),
        base_template='admin/base.html'
    )

    # Add model views
    # admin.add_view(FilterView(db_session, name='Filters', category='Configuration'))

    admin.add_view(ConfigSetView(db_session, name='Configuration Sets', endpoint='config_sets', category='Configuration'))
    admin.add_view(StationView(db_session, name='Stations', category='Configuration'))
    admin.add_view(DataAvailabilityView(db_session, name='Data Availability', category='Data'))

    admin.add_view(WorkflowStepsView(db_session, name='Workflow Steps', category='Workflows'))

    admin.add_view(JobView(db_session, name='Jobs', category='Processing'))

    admin.add_view(ConfigView(db_session, name='Configuration (Expert)', category='Expert'))
    # Add API view
    # admin.add_view(WorkflowAPIView(name='Workflow API', category='API'))

    # Add custom templates
    # app.jinja_env.globals['WORKFLOW_CHAIN_BUILDER_TEMPLATE'] = WORKFLOW_CHAIN_BUILDER_TEMPLATE
    app.jinja_env.globals['CONFIG_BULK_UPDATE_TEMPLATE'] = CONFIG_BULK_UPDATE_TEMPLATE
    app.jinja_env.globals['MSNOISE_INDEX_TEMPLATE'] = MSNOISE_INDEX_TEMPLATE


    return admin


# ============================================================================
# Main Entry Point
# ============================================================================

def main(port=5099):
    """Main entry point for the MSNoise admin interface"""

    logger.info("Starting MSNoise Admin Interface")

    # Create admin interface
    admin = create_admin_app()

    # Configure Flask app
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    # Add template filters
    @app.template_filter('datetime')
    def datetime_filter(value):
        if value is None:
            return ''
        return value.strftime('%Y-%m-%d %H:%M:%S')

    # Error handlers
    @app.errorhandler(404)
    def not_found(error):
        return render_template_string('<h1>Page not found</h1>'), 404

    @app.errorhandler(500)
    def internal_error(error):
        db_session.rollback()
        return render_template_string('<h1>Internal server error</h1>'), 500

    # Cleanup
    @app.teardown_appcontext
    def shutdown_session(exception=None):
        db_session.close()

    # Start the application
    try:
        app.run(host='0.0.0.0', port=port, debug=True)
    except KeyboardInterrupt:
        logger.info("Shutting down MSNoise Admin Interface")
    finally:
        db_session.close()


if __name__ == '__main__':
    main()
